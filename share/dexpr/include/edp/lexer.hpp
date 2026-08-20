#ifndef EDP_LEXER_HPP
#define EDP_LEXER_HPP

#include <edp/tokens.hpp>
#include <string>

namespace edp {

class Lexer {

public:
  explicit Lexer(std::string input);
  ~Lexer() = default;

  Token next_token();

private:
  std::string input_;
  int position_;
  int read_position_;
  char current_char_;

  // functions
  void skip_whitespace();

  std::string read_identifier();
  // seen_dot is true when the caller already consumed a leading '.', so that
  // a second one is not folded into the same literal.
  std::string read_number(bool seen_dot = false);
  bool read_to_delim(char ch, std::string& out);

  char peek_char() const;
  void read_char();
  Token make_token(TokenTypes kind) const;

};

}

#endif
