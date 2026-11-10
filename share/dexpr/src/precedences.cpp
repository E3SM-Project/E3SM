#include <dexpr/precedences.hpp>
#include <dexpr/tokens.hpp>
#include <stdexcept>

namespace dexpr::parser {

// Lowest,
// Equal,
// LessGreater,
// Sum,
// Product,
// Prefix,
// Bounds,
// Call,

Precedence token_precedence(TokenTypes type) {
  switch (type) {

  case TokenTypes::Assign:
    return Precedence::Assignment;

  case TokenTypes::And:
  case TokenTypes::Or:
    return Precedence::Logical;

  case TokenTypes::Equal:
  case TokenTypes::NotEqual:
    return Precedence::Equalitative;

  case TokenTypes::GreaterThan:
  case TokenTypes::GreaterEqual:
  case TokenTypes::LessThan:
  case TokenTypes::LessEq:
    return Precedence::Comparison;

  case TokenTypes::Plus:
  case TokenTypes::Minus:
    return Precedence::Additive;

  case TokenTypes::Slash:
  case TokenTypes::Asterisk:
    return Precedence::Multiplicative;

  case TokenTypes::Bang:
    return Precedence::Prefix;

  case TokenTypes::Exp:
    return Precedence::Exponent;

  case TokenTypes::Dot:
  case TokenTypes::LeftParen:
    return Precedence::Call;

  default:
    return Precedence::Lowest;
  }
}

Precedence right_binding_precedence(TokenTypes type) {
  const auto prec = token_precedence(type);

  switch (type) {
  // '**' is the only right-associative operator: 2**3**2 is 2**(3**2).
  case TokenTypes::Exp:
    return static_cast<Precedence>(static_cast<int>(prec) - 1);
  default:
    return prec;
  }
}

} // namespace dexpr::parser
