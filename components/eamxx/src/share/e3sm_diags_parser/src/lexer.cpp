#include <edp/lexer.hpp>
#include <edp/tokens.hpp>
#include <algorithm>
#include <cctype>
#include <utility>

namespace {

bool is_valid_identifier(const char ch) {
  return std::isalpha(static_cast<unsigned char>(ch)) || ch == '_';
}

bool is_numeric(const char ch) {
  return std::isdigit(static_cast<unsigned char>(ch));
}

} // namespace

namespace edp {

Lexer::Lexer(std::string input)
    : input_{std::move(input)}, position_{0}, read_position_{0},
      current_char_{'\0'} {
  // The input is deliberately NOT case-folded here. Folding the whole buffer
  // made keywords case-insensitive, but it also rewrote string literals --
  // which are data, not syntax -- so 'MyVar' silently became 'myvar'. Case
  // insensitivity is applied where it belongs instead: identifier_lookup()
  // folds before matching keywords, and read_number() accepts either 'e' or
  // 'E' as the exponent marker.
  read_char();
}

void Lexer::read_char() {
  if (read_position_ >= input_.length()) {
    current_char_ = '\0';
  } else {
    current_char_ = input_.at(read_position_);
  }
  position_ = read_position_;
  read_position_ += 1;
}

char Lexer::peek_char() const {
  if (read_position_ >= input_.length()) {
    return '\0';
  } else {
    return input_[read_position_];
  }
}

void Lexer::skip_whitespace() {
  while (std::isspace(static_cast<unsigned char>(current_char_))) {
    read_char();
  }
}

Token Lexer::make_token(TokenTypes kind) const {
  return {kind, std::string(1, current_char_)};
}

bool Lexer::read_to_delim(char ch, std::string& out) {
  const auto start_pos = position_+1;
  while (peek_char() != ch && peek_char() != '\0') {
    read_char();
  }

  const bool closed = peek_char() == ch;
  read_char();

  auto count = position_ - start_pos;
  out = input_.substr(start_pos, count);
  return closed;
}

std::string Lexer::read_number() {
  const auto start_pos = position_;

  // At most one decimal point belongs to a single number: "1.2.3" is two
  // numbers, not one malformed one
  bool seen_dot = false;
  while (is_numeric(current_char_) || (current_char_ == '.' && !seen_dot)) {
    if (current_char_ == '.') {
      seen_dot = true;
    }
    read_char();
  }

  // Only consume an exponent if a digit actually follows it
  if (current_char_ == 'e' || current_char_ == 'E') {
    const auto sign_offset =
        (peek_char() == '+' || peek_char() == '-') ? 1 : 0;
    const auto digit_pos = read_position_ + sign_offset;

    if (digit_pos < static_cast<int>(input_.length()) &&
        is_numeric(input_[digit_pos])) {
      read_char(); // 'e'
      if (current_char_ == '+' || current_char_ == '-') {
        read_char(); // sign
      }
      while (is_numeric(current_char_)) {
        read_char();
      }
    }
  }

  return input_.substr(start_pos, position_ - start_pos);
}

std::string Lexer::read_identifier() {
  auto start_pos = position_;
  while (is_valid_identifier(current_char_)) {
    read_char();
  }
  auto length = position_ - start_pos;
  return input_.substr(start_pos, length);
}

Token Lexer::next_token() {

  skip_whitespace();

  Token tok;

  switch (current_char_) {
  case '=':
    if (peek_char() == '=') {
      tok = {TokenTypes::Equal, "=="};
      read_char();
    } else {
      tok = make_token(TokenTypes::Assign);
    }
    break;
  case '(':
    tok = make_token(TokenTypes::LeftParen);
    break;
  case ')':
    tok = make_token(TokenTypes::RightParen);
    break;
  case ',':
    tok = make_token(TokenTypes::Comma);
    break;
  case '+':
    tok = make_token(TokenTypes::Plus);
    break;
  case '-':
    tok = make_token(TokenTypes::Minus);
    break;
  case '[':
    tok = make_token(TokenTypes::ArrayLeftBracket);
    break;
  case ']':
    tok = make_token(TokenTypes::ArrayRightBracket);
    break;
  case '/':
    tok = make_token(TokenTypes::Slash);
    break;
  case '*':
    if (peek_char() == '*') {
      read_char();
      tok = {TokenTypes::Exp, "**"};
    } else {
      tok = make_token(TokenTypes::Asterisk);
    }
    break;
  case '<':
    if (peek_char() == '=') {
      read_char();
      tok = {TokenTypes::LessEq, "<="};
    } else {
      tok = make_token(TokenTypes::LessThan);
    }
    break;

  case '>':
    if (peek_char() == '=') {
      read_char();
      tok = {TokenTypes::GreaterEqual, ">="};
    } else {
      tok = make_token(TokenTypes::GreaterThan);
    }
    break;
  case '!':
    if (peek_char() == '=') {
      read_char();
      tok = {TokenTypes::NotEqual, "!="};
    } else {
      // '!' is an alias for the `not` keyword
      tok = make_token(TokenTypes::Bang);
    }
    break;
  case '\0':
    return {TokenTypes::EndofFile, ""};
  case ':':
    tok = make_token(TokenTypes::Colon);
    break;
  case '"':
  case '\'': {
    std::string literal;
    const bool closed = read_to_delim(current_char_, literal);
    // An unterminated literal is reported as Illegal rather than silently
    // accepted as a well-formed string.
    tok = {closed ? TokenTypes::String : TokenTypes::Illegal, literal};
    break;
  }
  case '.': {
    if (is_numeric(peek_char())) {
      read_char();
      auto number = read_number();
      number.insert(0, "0.");
      return {TokenTypes::Float, number};
    } else {
        tok = make_token(TokenTypes::Dot);
        break;
      }
  }
  default: {

    if (is_valid_identifier(current_char_)) {
      return identifier_lookup({TokenTypes::Identifier, read_identifier()});
    } else if (is_numeric(current_char_)) {
      auto number = read_number();

      if (number.find_first_of(".e") != std::string::npos) {
        return {TokenTypes::Float, number};
      } else {
        return {TokenTypes::Integer, number};
      }
    } else {
      auto illegal = make_token(TokenTypes::Illegal);
      read_char();
      return illegal;
    }

  } // default
  } // end switch
  read_char();
  return tok;
}

} // namespace edp
