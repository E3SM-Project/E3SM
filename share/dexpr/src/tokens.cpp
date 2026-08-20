#include <edp/tokens.hpp>
#include <algorithm>
#include <cctype>
#include <stdexcept>
#include <string>

namespace edp {

std::string_view to_string(TokenTypes type) {
  switch (type) {

  case TokenTypes::EndofFile:
    return "EndofFile";
  case TokenTypes::Illegal:
    return "Illegal";

  case TokenTypes::Identifier:
    return "Identifer";
  case TokenTypes::Integer:
    return "Integer";
  case TokenTypes::Float:
    return "Float";
  case TokenTypes::String:
    return "String";
  // Operators
  case TokenTypes::Assign:
    return "Assign";
  case TokenTypes::Plus:
    return "Plus";
  case TokenTypes::Minus:
    return "Minus";
  case TokenTypes::Asterisk:
    return "Asterisk";
  case TokenTypes::Bang:
    return "Bang";
  case TokenTypes::Slash:
    return "Slash";
  case TokenTypes::Exp:
    return "Exp";
  case TokenTypes::Equal:
    return "Equal";
  case TokenTypes::GreaterThan:
    return "GreaterThan";
  case TokenTypes::GreaterEqual:
    return "GreaterEqual";
  case TokenTypes::LessThan:
    return "LessThan";
  case TokenTypes::LessEq:
    return "LessEq";
  case TokenTypes::NotEqual:
    return "NotEqual";
  case TokenTypes::Or:
    return "Or";
  case TokenTypes::And:
    return "And";
  case TokenTypes::Concat:
    return "Concat";
  case TokenTypes::Dot:
    return "Dot";

  // DELIMITERS
  case TokenTypes::Comma:
    return "Comma";
  case TokenTypes::LeftParen:
    return "LeftParen";
  case TokenTypes::RightParen:
    return "RightParen";
  case TokenTypes::Colon:
    return "Colon";
  case TokenTypes::Semicolon:
    return "Semicolon";
  case TokenTypes::Percent:
    return "Percent";
  case TokenTypes::DoubleColon:
    return "DoubleColon";
  case TokenTypes::ArrayLeftBracket:
    return "ArrayLeftBracket";
  case TokenTypes::ArrayRightBracket:
    return "ArrayRightBracket";
  default:
    return "UNKNOWN";
  }
}

std::string to_string(const Token& tok) {
  return "{Type: " + std::string(to_string(tok.type)) +
         ", Literal: " + tok.literal + "}";
}
Token identifier_lookup(const Token& tok) {
  // This function checks to see if an identifier is a keyword
  // keywords are case-insensitive; identifiers are case-sensitive
  // identifier case intact.
  std::string folded = tok.literal;
  std::transform(
      folded.begin(), folded.end(), folded.begin(),
      [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

  if (const auto it = keywords.find(folded); it != keywords.end()) {
    return it->second;
  }
  return tok;
}

std::ostream& operator<<(std::ostream& os, const Token& tok) {
  os << to_string(tok);
  return os;
}
std::string binary_op_to_string(const TokenTypes type) {

  switch (type) {
  case TokenTypes::Plus:
    return "+";
  case TokenTypes::Minus:
    return "-";
  case TokenTypes::Asterisk:
    return "*";
  case TokenTypes::Slash:
    return "/";
  case TokenTypes::Equal:
    return "==";
  case TokenTypes::NotEqual:
    return "!=";
  case TokenTypes::LessThan:
    return "<";
  case TokenTypes::GreaterThan:
    return ">";
  case TokenTypes::Assign:
    return "=";
  case TokenTypes::Exp:
    return "**";
  case TokenTypes::Dot:
    return ".";
  default:
    throw std::invalid_argument{"Invalid Binary Operator" + std::string(to_string(type))};
  }
}

std::string unary_op_to_string(const TokenTypes type) {
  switch (type) {
  case TokenTypes::Plus:
    return "+";
  case TokenTypes::Minus:
    return "-";
  case TokenTypes::Bang:
    return "!";
  default:
    throw std::invalid_argument{"Invalid Unary Operator"};
  }
}
}
