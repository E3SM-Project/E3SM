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

  std::transform(
      input_.begin(), input_.end(), input_.begin(),
      [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
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

std::string Lexer::read_to_delim(char ch) {
  auto start_pos = position_+1;
  while (peek_char() != ch) {
    read_char();
  }
  read_char();
  auto count = position_ - start_pos;
  return input_.substr(start_pos, count);
}

std::string Lexer::read_number() {
  auto start_pos = position_;
  while (is_numeric(current_char_)) {
    read_char();
    if (current_char_ == '.') {
      read_char();
    }
  }
  check_precision();
  return input_.substr(start_pos, position_ - start_pos);
}

void Lexer::check_precision() {
  const auto peek_ch = peek_char();
  if (std::isspace(static_cast<unsigned char>(peek_ch))) {
    return;
  }
  if (current_char_ == 'e') {
    read_char();
    if (current_char_ == '+' || current_char_ == '-') {
      read_char();
    }
    read_char();
    read_number();
  }
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
  case '\n':
    tok = make_token(TokenTypes::Newline);
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
  case '\0':
    return {TokenTypes::EndofFile, ""};
  case ':':
    tok = make_token(TokenTypes::Colon);
    break;
  case '"':
  case '\'':
    tok = {TokenTypes::String, read_to_delim(current_char_)};
    break;
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

      if (number.find('e') != std::string::npos) {
        return {TokenTypes::Float, number};
      } else {
        return {TokenTypes::Integer, number};
      }
    } else {
      return make_token(TokenTypes::Illegal);
    }

  } // default
  } // end switch
  read_char();
  return tok;
}

} // namespace edp
