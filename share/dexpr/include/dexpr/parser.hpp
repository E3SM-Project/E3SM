#ifndef DEXPR_PARSER_HPP
#define DEXPR_PARSER_HPP

#include <dexpr/ast.hpp>
#include <dexpr/lexer.hpp>
#include <dexpr/precedences.hpp>
#include <dexpr/tokens.hpp>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace dexpr::parser {

class ParserError : public std::runtime_error {
public:
  explicit ParserError(const std::vector<std::string>& errors)
      : std::runtime_error(join_msgs(errors)) {}

private:
  static std::string join_msgs(const std::vector<std::string>& errors) {
    std::string result = "Parser errors:\n";

    for (const auto& error : errors) {
      result += "  - " + error + '\n';
    }

    return result;
  }
};

class Parser {

public:
  explicit Parser(Lexer lexer);

  ast::ExprPtr parse();
  bool has_errors();

private:
  Lexer lexer_;
  Token cur_token_;
  Token peek_token_;
  std::vector<std::string> errors_;

  // Prefix parsing functions must return ExprPtr and take no arugments (other
  // than this)
  using PrefixFn = ast::ExprPtr (Parser::*)();
  // Infix parsing function must return ExprPtr and take an ExprPtr as an
  // argument
  using InfixFn = ast::ExprPtr (Parser::*)(ast::ExprPtr);

  std::unordered_map<TokenTypes, PrefixFn> prefix_parse_fns_;
  std::unordered_map<TokenTypes, InfixFn> infix_parse_fns_;

  // Functions
  void add_error(std::string msg);
  void next_token();

  bool peek_token_is(TokenTypes expected_type);

  bool expect_peek_and_advance(TokenTypes expected_type);
  Precedence peek_precedence();

  ast::ExprPtr parse_expression(Precedence prec);

  // prefix member functions:
  ast::ExprPtr parse_identifier();
  ast::ExprPtr parse_integer_literal();
  ast::ExprPtr parse_string_literal();
  ast::ExprPtr parse_float_literal();
  ast::ExprPtr parse_prefix_expression();
  ast::ExprPtr parse_grouped_expression();
  ast::ExprPtr parse_array_expression();

  // infix member functions:
  ast::ExprPtr parse_infix_expression(ast::ExprPtr left_expr);
  ast::ExprPtr parse_function_expression(ast::ExprPtr expr);
  std::vector<ast::ExprPtr> parse_list_of_expressions(TokenTypes end_token);
};

} // namespace dexpr::parser

#endif
