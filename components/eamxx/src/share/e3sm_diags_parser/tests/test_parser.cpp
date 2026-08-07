#include <edp/ast.hpp>
#include "catch2/catch_message.hpp"
#include <edp/lexer.hpp>
#include <edp/parser.hpp>
#include <edp/tokens.hpp>
#include <catch2/catch_test_macros.hpp>
#include <iostream>

namespace edp {

TEST_CASE("Test Parse expressions") {
  std::string input = "x*y.derivative(dx=dy,['col']).where(x>0)";

  parser::Parser parser{Lexer{input}};

  auto expr = parser.parse();
  auto str_ = to_string(*expr);
  std::cout << "Parsed Expression: \n" << str_;
  INFO("Parsed Expression: \n" << str_);
  CHECK(str_ == "(x*((y.derivative((dx=dy), ['col'])).where((x>0))))");
}

} // namespace edp
