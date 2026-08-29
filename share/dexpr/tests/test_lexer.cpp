#include <catch2/catch_test_macros.hpp>

#include <dexpr/lexer.hpp>
#include <dexpr/tokens.hpp>

namespace dexpr
{

namespace
{ // anonymous

// safety net inside tests
constexpr std::size_t k_max_token = 1024;

std::vector<Token>
lex_all(std::string input)
{
  Lexer lexer{input};
  std::vector<Token> tokens;
  while (tokens.size() < k_max_token) {
    auto tok = lexer.next_token();
    tokens.push_back(tok);
    if (tok.type == TokenTypes::EndofFile) {
      break;
    }
  }
  return tokens;
}

void
check_tokens(const std::string &input, const std::vector<Token> &expected)
{
  const auto actual = lex_all(input);

  INFO("Input: " << input);
  REQUIRE(actual.size() < k_max_token);
  REQUIRE(actual.size() == expected.size());

  for (std::size_t i = 0; i < expected.size(); ++i) {
    INFO("Token #" << i << "\n  Expected: " << to_string(expected[i])
                   << "\n  Received: " << to_string(actual[i]));
    REQUIRE(expected[i].type == actual[i].type);     // unwind
    CHECK(expected[i].literal == actual[i].literal); // keep going
  }
}

const Token k_eof{TokenTypes::EndofFile, ""};

} // namespace

TEST_CASE("lexer: empty and whitespace-only input", "[lexer]")
{
  check_tokens("", {k_eof});
  check_tokens("   \t  ", {k_eof});
}

// TODO: for now, newlines are simply ignored as whitespace; tacky but workable
TEST_CASE("lexer: newline is like whitespace", "[lexer]")
{
  check_tokens("x\ny _w my_z\n", {
                                     {TokenTypes::Identifier, "x"},
                                     {TokenTypes::Identifier, "y"},
                                     {TokenTypes::Identifier, "_w"},
                                     {TokenTypes::Identifier, "my_z"},
                                     k_eof,
                                 });
}

TEST_CASE("lexer: single-character operators and delimiters", "[lexer]")
{
  check_tokens("+-*/(),:[].", {
                                  {TokenTypes::Plus, "+"},
                                  {TokenTypes::Minus, "-"},
                                  {TokenTypes::Asterisk, "*"},
                                  {TokenTypes::Slash, "/"},
                                  {TokenTypes::LeftParen, "("},
                                  {TokenTypes::RightParen, ")"},
                                  {TokenTypes::Comma, ","},
                                  {TokenTypes::Colon, ":"},
                                  {TokenTypes::ArrayLeftBracket, "["},
                                  {TokenTypes::ArrayRightBracket, "]"},
                                  {TokenTypes::Dot, "."},
                                  k_eof,
                              });
}

TEST_CASE("lexer: multi-character operators", "[lexer]")
{
  check_tokens("== <= >= ** < > =", {
                                        {TokenTypes::Equal, "=="},
                                        {TokenTypes::LessEq, "<="},
                                        {TokenTypes::GreaterEqual, ">="},
                                        {TokenTypes::Exp, "**"},
                                        {TokenTypes::LessThan, "<"},
                                        {TokenTypes::GreaterThan, ">"},
                                        {TokenTypes::Assign, "="},
                                        k_eof,
                                    });
}

TEST_CASE("lexer: integer literals", "[lexer]")
{
  check_tokens("0 5 42 0 55", {
                                  {TokenTypes::Integer, "0"},
                                  {TokenTypes::Integer, "5"},
                                  {TokenTypes::Integer, "42"},
                                  {TokenTypes::Integer, "0"},
                                  {TokenTypes::Integer, "55"},
                                  k_eof,
                              });
}

TEST_CASE("lexer: decimal float literals", "[lexer]")
{
  check_tokens("1.", {{TokenTypes::Float, "1."}, k_eof});
  check_tokens("1.5", {{TokenTypes::Float, "1.5"}, k_eof});
  check_tokens("0.25", {{TokenTypes::Float, "0.25"}, k_eof});
  check_tokens(".025", {{TokenTypes::Float, "0.025"}, k_eof});
}

TEST_CASE("lexer: exponent-form float literals", "[lexer]")
{
  check_tokens("1.0e-4", {{TokenTypes::Float, "1.0e-4"}, k_eof});
  check_tokens("2.5e+3", {{TokenTypes::Float, "2.5e+3"}, k_eof});
  check_tokens("1.5e3", {{TokenTypes::Float, "1.5e3"}, k_eof});
  check_tokens("1e5", {{TokenTypes::Float, "1e5"}, k_eof});
  check_tokens("1.e5", {{TokenTypes::Float, "1.e5"}, k_eof});
  // Either case marks an exponent, and the literal keeps whichever was
  // written -- the lexer no longer rewrites the input. std::stof accepts both.
  check_tokens("1.0E-4", {{TokenTypes::Float, "1.0E-4"}, k_eof});
}

TEST_CASE("lexer: a number followed by an operator", "[lexer]")
{
  check_tokens("1.0e-4+2", {
                               {TokenTypes::Float, "1.0e-4"},
                               {TokenTypes::Plus, "+"},
                               {TokenTypes::Integer, "2"},
                               k_eof,
                           });
}

TEST_CASE("lexer: truncated exponent does not throw", "[lexer]")
{
  check_tokens("1.0e", {
                           {TokenTypes::Float, "1.0"},
                           {TokenTypes::Identifier, "e"},
                           k_eof,
                       });
  check_tokens("1.0e+", {
                            {TokenTypes::Float, "1.0"},
                            {TokenTypes::Identifier, "e"},
                            {TokenTypes::Plus, "+"},
                            k_eof,
                        });
}

TEST_CASE("lexer: a second decimal point ends the number", "[lexer]")
{
  check_tokens("1.2.3", {
                            {TokenTypes::Float, "1.2"},
                            {TokenTypes::Float, "0.3"},
                            k_eof,
                        });
}

TEST_CASE("lexer: identifiers may contain digits", "[lexer]")
{
  check_tokens("bc_a1 so4_a2 O3 num_a1", {
                                             {TokenTypes::Identifier, "bc_a1"},
                                             {TokenTypes::Identifier, "so4_a2"},
                                             {TokenTypes::Identifier, "O3"},
                                             {TokenTypes::Identifier, "num_a1"},
                                             k_eof,
                                         });
  check_tokens("_a1", {{TokenTypes::Identifier, "_a1"}, k_eof});
  check_tokens("dst_a3_at_lev_10", {{TokenTypes::Identifier, "dst_a3_at_lev_10"}, k_eof});
}

TEST_CASE("lexer: an identifier cannot start with a digit", "[lexer]")
{
  check_tokens("2x", {
                         {TokenTypes::Integer, "2"},
                         {TokenTypes::Identifier, "x"},
                         k_eof,
                     });
  check_tokens("500hPa", {
                             {TokenTypes::Integer, "500"},
                             {TokenTypes::Identifier, "hPa"},
                             k_eof,
                         });
}

TEST_CASE("lexer: digits do not break member access", "[lexer]")
{
  check_tokens("bc_a1.mean", {
                                 {TokenTypes::Identifier, "bc_a1"},
                                 {TokenTypes::Dot, "."},
                                 {TokenTypes::Identifier, "mean"},
                                 k_eof,
                             });
}

TEST_CASE("lexer: a keyword with a digit appended is an identifier", "[lexer]")
{
  check_tokens("and2 or1 not3", {
                                    {TokenTypes::Identifier, "and2"},
                                    {TokenTypes::Identifier, "or1"},
                                    {TokenTypes::Identifier, "not3"},
                                    k_eof,
                                });
}

TEST_CASE("lexer: keywords are case-insensitive, normalized but not identifiers", "[lexer]")
{
  check_tokens("AND or Not  NOT android nothing", {
                                                      {TokenTypes::And, "and"},
                                                      {TokenTypes::Or, "or"},
                                                      {TokenTypes::Bang, "!"},
                                                      {TokenTypes::Bang, "!"},
                                                      {TokenTypes::Identifier, "android"},
                                                      {TokenTypes::Identifier, "nothing"},
                                                      k_eof,
                                                  });
}

TEST_CASE("lexer: string literals in both quote styles", "[lexer]")
{
  check_tokens("'col'", {{TokenTypes::String, "col"}, k_eof});
  check_tokens("\"col\"", {{TokenTypes::String, "col"}, k_eof});
}

TEST_CASE("lexer: string literals preserve case", "[lexer]")
{
  check_tokens("'MyVar'", {{TokenTypes::String, "MyVar"}, k_eof});
  check_tokens("\"SHOC_tke\"", {{TokenTypes::String, "SHOC_tke"}, k_eof});
  check_tokens("MyField T_mid", {
                                    {TokenTypes::Identifier, "MyField"},
                                    {TokenTypes::Identifier, "T_mid"},
                                    k_eof,
                                });
}

TEST_CASE("lexer: unterminated string literal is illegal", "[lexer]")
{
  check_tokens("'abc", {{TokenTypes::Illegal, "abc"}, k_eof});
  check_tokens("\"abc", {{TokenTypes::Illegal, "abc"}, k_eof});
  check_tokens("'", {{TokenTypes::Illegal, ""}, k_eof});
  check_tokens("'ok' + 'bad", {
                                  {TokenTypes::String, "ok"},
                                  {TokenTypes::Plus, "+"},
                                  {TokenTypes::Illegal, "bad"},
                                  k_eof,
                              });
}

// The Illegal branch in next_token() returned early, skipping the read_char()
// at the end of the function -- so the scanner stayed parked on the offending
// character and re-emitted it forever. Any consumer draining to EndofFile
// (lex_all included) hung rather than failing.
TEST_CASE("lexer: illegal characters advance the scanner", "[lexer]")
{
  check_tokens("a@b", {
                          {TokenTypes::Identifier, "a"},
                          {TokenTypes::Illegal, "@"},
                          {TokenTypes::Identifier, "b"},
                          k_eof,
                      });
  check_tokens("@", {{TokenTypes::Illegal, "@"}, k_eof});
  check_tokens("@#@", {
                          {TokenTypes::Illegal, "@"},
                          {TokenTypes::Illegal, "#"},
                          {TokenTypes::Illegal, "@"},
                          k_eof,
                      });
}

TEST_CASE("lexer: not-equal operator", "[lexer]")
{
  check_tokens("3 != 4", {
                             {TokenTypes::Integer, "3"},
                             {TokenTypes::NotEqual, "!="},
                             {TokenTypes::Integer, "4"},
                             k_eof,
                         });
  check_tokens("x!=y", {
                           {TokenTypes::Identifier, "x"},
                           {TokenTypes::NotEqual, "!="},
                           {TokenTypes::Identifier, "y"},
                           k_eof,
                       });
}

TEST_CASE("lexer: bang is an alias for not", "[lexer]")
{
  check_tokens("!x", {
                         {TokenTypes::Bang, "!"},
                         {TokenTypes::Identifier, "x"},
                         k_eof,
                     });
  check_tokens("not x", {
                            {TokenTypes::Bang, "!"},
                            {TokenTypes::Identifier, "x"},
                            k_eof,
                        });
}

TEST_CASE("lexer: bang alias does not swallow not-equal", "[lexer]")
{
  check_tokens("x != y", {
                             {TokenTypes::Identifier, "x"},
                             {TokenTypes::NotEqual, "!="},
                             {TokenTypes::Identifier, "y"},
                             k_eof,
                         });
  check_tokens("a ! = b", {
                              {TokenTypes::Identifier, "a"},
                              {TokenTypes::Bang, "!"},
                              {TokenTypes::Assign, "="},
                              {TokenTypes::Identifier, "b"},
                              k_eof,
                          });
  check_tokens("! x", {
                          {TokenTypes::Bang, "!"},
                          {TokenTypes::Identifier, "x"},
                          k_eof,
                      });
}

// A malformed literal must not be quietly repaired into a valid-looking one:
// ".5.3" used to lex as a single Float("0.5.3"), which the parser then read as
// 0.5 and silently dropped the rest.
TEST_CASE("lexer: a second decimal point starts a new number", "[lexer]")
{
  check_tokens(".5.3", {
                           {TokenTypes::Float, "0.5"},
                           {TokenTypes::Float, "0.3"},
                           k_eof,
                       });
  // The leading-dot form must agree with the ordinary form
  check_tokens("1.2.3", {
                            {TokenTypes::Float, "1.2"},
                            {TokenTypes::Float, "0.3"},
                            k_eof,
                        });
  check_tokens(".5", {{TokenTypes::Float, "0.5"}, k_eof});
}

// "1E5" used to lex as Integer("1E5"); integer parsing then stopped at the 'E'
// and the expression silently evaluated to 1.
TEST_CASE("lexer: exponent marker is case-insensitive", "[lexer]")
{
  check_tokens("1E5", {{TokenTypes::Float, "1E5"}, k_eof});
  check_tokens("1e5", {{TokenTypes::Float, "1e5"}, k_eof});
  check_tokens("1E+5", {{TokenTypes::Float, "1E+5"}, k_eof});
  check_tokens("1E-5", {{TokenTypes::Float, "1E-5"}, k_eof});
  check_tokens("2.5E3", {{TokenTypes::Float, "2.5E3"}, k_eof});
  // Still an identifier when no digit follows, not a broken number
  check_tokens("1E", {
                         {TokenTypes::Integer, "1"},
                         {TokenTypes::Identifier, "E"},
                         k_eof,
                     });
}

TEST_CASE("lexer: tokens carry their source position", "[lexer]")
{
  //                          1234567
  const auto tokens = lex_all("x + 12");
  REQUIRE(tokens.size() == 4);
  CHECK(tokens[0].line == 1);
  CHECK(tokens[0].column == 1); // x
  CHECK(tokens[1].column == 3); // +
  CHECK(tokens[2].column == 5); // 12
  CHECK(tokens[3].column == 7); // end of input
}

TEST_CASE("lexer: position survives multi-character tokens", "[lexer]")
{
  //                          123456789
  const auto tokens = lex_all("ab <= 'c'");
  REQUIRE(tokens.size() == 4);
  CHECK(tokens[0].column == 1); // ab
  CHECK(tokens[1].column == 4); // <=
  CHECK(tokens[2].column == 7); // 'c'
}

TEST_CASE("lexer: newlines advance the line and reset the column", "[lexer]")
{
  const auto tokens = lex_all("a +\n  b\nc");
  REQUIRE(tokens.size() == 5);
  CHECK(tokens[0].line == 1);
  CHECK(tokens[0].column == 1); // a
  CHECK(tokens[1].line == 1);
  CHECK(tokens[1].column == 3); // +
  CHECK(tokens[2].line == 2);
  CHECK(tokens[2].column == 3); // b, past two spaces
  CHECK(tokens[3].line == 3);
  CHECK(tokens[3].column == 1); // c
}

TEST_CASE("some lexer token stream", "[lexer]")
{
  check_tokens(" not x <= 1.0e-4 and y + 5=1", {
                                                   {TokenTypes::Bang, "!"},
                                                   {TokenTypes::Identifier, "x"},
                                                   {TokenTypes::LessEq, "<="},
                                                   {TokenTypes::Float, "1.0e-4"},
                                                   {TokenTypes::And, "and"},
                                                   {TokenTypes::Identifier, "y"},
                                                   {TokenTypes::Plus, "+"},
                                                   {TokenTypes::Integer, "5"},
                                                   {TokenTypes::Assign, "="},
                                                   {TokenTypes::Integer, "1"},
                                                   {TokenTypes::EndofFile, ""},
                                               });
}

} // namespace dexpr
