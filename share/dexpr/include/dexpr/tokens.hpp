/**
 * @file tokens.hpp
 * @brief Defines the all the tokens that may be produced by the lexer and
 * consumed by the Parser
 */
#ifndef DEXPR_TOKENS_HPP
#define DEXPR_TOKENS_HPP

#include <array>
#include <ostream>
#include <string>
#include <string_view>
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
  Dot,

  // DELIMITERS
  Comma,
  LeftParen,
  RightParen,
  Colon,
  ArrayLeftBracket,
  ArrayRightBracket,
};

std::string_view to_string(TokenTypes type);

struct Token {
  TokenTypes type;
  std::string literal;
  // Start of the token in the input, 1-based.
  int line = 1;
  int column = 1;
};

// "line 1, column 7"
std::string position_of(const Token& tok);
std::string to_string(const Token& tok);

// constexpr rather than a namespace-scope const container: the latter has
// internal linkage, so it would be built once per translation unit.
struct Keyword {
  std::string_view name; // spelling, always lower case; matching folds first
  TokenTypes type;
  std::string_view literal; // what the resulting token carries
};

inline constexpr std::array<Keyword, 3> keywords{{
    {"or", TokenTypes::Or, "or"},
    {"and", TokenTypes::And, "and"},
    {"not", TokenTypes::Bang, "!"},
}};
Token identifier_lookup(const Token& tok);

std::string binary_op_to_string(const TokenTypes type);
std::string unary_op_to_string(const TokenTypes type);

} // namespace dexpr

#endif
