#include <dexpr/tokens.hpp>
#include <dexpr/precedences.hpp>
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

  case TokenTypes::Equal:
  case TokenTypes::NotEqual:
  case TokenTypes::Assign:
  case TokenTypes::And:
  case TokenTypes::Or:
    return Precedence::Equal;

  case TokenTypes::GreaterThan:
  case TokenTypes::GreaterEqual:
  case TokenTypes::LessThan:
  case TokenTypes::LessEq:
    return Precedence::LessGreater;

  case TokenTypes::Plus:
  case TokenTypes::Minus:
    return Precedence::Sum;

  case TokenTypes::Slash:
  case TokenTypes::Asterisk:
    return Precedence::Product;

  case TokenTypes::Bang:
    return Precedence::Prefix;

  case TokenTypes::Exp:
    return Precedence::Exponent;

  case TokenTypes::Colon:
    return Precedence::Bounds;

  case TokenTypes::Dot:
  case TokenTypes::LeftParen:
    return Precedence::Call;

  default:
    return Precedence::Lowest;
  }

}

Precedence cur_precedence(TokenTypes type) {
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
