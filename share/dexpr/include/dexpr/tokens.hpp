#ifndef DEXPR_TOKENS_HPP
#define DEXPR_TOKENS_HPP

#include <ostream>
#include <string>
#include <string_view>
#include <unordered_map>
namespace dexpr {

enum class TokenTypes {
  EndofFile,
  Illegal,

  Identifier,
  Integer,
  Float,
  String,

  // Operators
  Assign,
  Plus,
  Minus,
  Asterisk,
  Bang,
  Slash,
  Exp,
  Equal,
  GreaterThan,
  GreaterEqual,
  LessThan,
  LessEq,
  NotEqual,
  Or,
  And,
  Concat,
  Dot,

  // DELIMITERS
  Comma,
  LeftParen,
  RightParen,
  Colon,
  Semicolon,
  Percent,
  DoubleColon,
  ArrayLeftBracket,
  ArrayRightBracket,
};

std::string_view to_string(TokenTypes type);

struct Token {
  TokenTypes type;
  std::string literal;
};
std::string to_string(const Token& tok);

const std::unordered_map<std::string, Token> keywords{
    {"or", {TokenTypes::Or, "or"}},
    {"and", {TokenTypes::And, "and"}},
    {"not", {TokenTypes::Bang, "!"}},
};
Token identifier_lookup(const Token& tok);

std::ostream& operator<<(std::ostream& os, const Token& tok);
std::string binary_op_to_string(const TokenTypes type);
std::string unary_op_to_string(const TokenTypes type) ;

} // namespace dexpr

#endif
