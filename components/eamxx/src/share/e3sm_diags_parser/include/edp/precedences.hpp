#ifndef EDP_PRECEDENCES_HPP
#define EDP_PRECEDENCES_HPP

#include <edp/tokens.hpp>
namespace edp::parser {

enum class Precedence {
  Lowest,
  Equal,
  LessGreater,
  Sum,
  Product,
  Prefix,
  Exponent,
  Bounds,
  Call,
};

Precedence token_precedence(TokenTypes type);
Precedence cur_precedence(TokenTypes type);
} // namespace edp::parser

#endif
