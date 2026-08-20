#ifndef DEXPR_PRECEDENCES_HPP
#define DEXPR_PRECEDENCES_HPP

#include <dexpr/tokens.hpp>
namespace dexpr::parser {

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
} // namespace dexpr::parser

#endif
