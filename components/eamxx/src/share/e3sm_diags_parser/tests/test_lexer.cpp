#include <catch2/catch_test_macros.hpp>

#include <edp/tokens.hpp>
#include <edp/lexer.hpp>

namespace edp {
TEST_CASE("lexer token stream") {
  Lexer lexer{" not x <= 1.0e-4 and y + 5=1"};

  const std::vector<Token> expected{
      {TokenTypes::Bang, "!"},     {TokenTypes::Identifier, "x"},
      {TokenTypes::LessEq, "<="},  {TokenTypes::Float, "1.0e-4"},
      {TokenTypes::And, "and"},    {TokenTypes::Identifier, "y"},
      {TokenTypes::Plus, "+"},     {TokenTypes::Integer, "5"},
      {TokenTypes::Assign, "="},   {TokenTypes::Integer, "1"},
      {TokenTypes::EndofFile, ""},
  };
  for (const auto& expected_token : expected) {
    auto my_token = lexer.next_token();
    INFO("Expected: " << to_string(expected_token)
                      << "\nReceived: " << to_string(my_token));
    CHECK(expected_token.type == my_token.type);
    CHECK(expected_token.literal == my_token.literal);
  }
}
} // namespace edp
