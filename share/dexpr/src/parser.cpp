#include <edp/parser.hpp>
#include <edp/ast.hpp>
#include <edp/precedences.hpp>
#include <edp/tokens.hpp>
#include <iostream>
#include <stdexcept>
#include <string>

namespace edp::parser {

bool Parser::cur_token_is(TokenTypes expected_type) {
  return cur_token_.type == expected_type;
};
bool Parser::peek_token_is(TokenTypes expected_type) {
  return peek_token_.type == expected_type;
};

bool Parser::expect_peek_and_advance(TokenTypes expected_type) {
  if (peek_token_is(expected_type)) {
    next_token();
    return true;
  } else {
    errors_.push_back("Expected " + std::string(to_string(expected_type)) +
                      "Got " + to_string(peek_token_));
    return false;
  }
}

Precedence Parser::cur_precedence() {
  return token_precedence(cur_token_.type);
}
Precedence Parser::peek_precedence() {
  return token_precedence(peek_token_.type);
}

bool Parser::has_errors() { return !errors_.empty(); }

void Parser::next_token() {
  cur_token_ = peek_token_;
  peek_token_ = lexer_.next_token();
  if (peek_token_is(TokenTypes::Illegal)) {
    std::cout << "Encountered Illegal Token: " << to_string(peek_token_);
  }
}

ast::ExprPtr Parser::parse_expression(Precedence prec) {
  const auto prefix = prefix_parse_fns_.find(cur_token_.type);
  if (prefix == prefix_parse_fns_.end()) {
    throw("Unexpected Prefix Token " + to_string(cur_token_));
  }
  const auto fn = prefix->second;
  auto left_expr = (this->*fn)();

  while (!peek_token_is(TokenTypes::EndofFile) && prec < peek_precedence()) {
    const auto infix_it = infix_parse_fns_.find(peek_token_.type);
    if (infix_it == infix_parse_fns_.end()) {
      return left_expr;
    }
    const auto infix_fn = infix_it->second;
    next_token();
    left_expr = (this->*infix_fn)(std::move(left_expr));
  }
  return left_expr;
}

ast::ExprPtr Parser::parse_identifier() {
  return ast::make_expression<ast::Identifier>(cur_token_.literal);
}
ast::ExprPtr Parser::parse_string_literal() {
  return ast::make_expression<ast::StringLiteral>(cur_token_.literal);
}
ast::ExprPtr Parser::parse_integer_literal() {
  return ast::make_expression<ast::IntegerLiteral>(
      std::stoi(cur_token_.literal));
}
ast::ExprPtr Parser::parse_float_literal() {
  return ast::make_expression<ast::FloatLiteral>(std::stof(cur_token_.literal));
}
ast::ExprPtr Parser::parse_prefix_expression() {
  auto op = cur_token_.type;
  auto right_expr = parse_expression(Precedence::Prefix);
  return ast::make_expression<ast::PrefixExpression>(op, std::move(right_expr));
}

ast::ExprPtr Parser::parse_grouped_expression() {
  next_token();
  auto expr = parse_expression(Precedence::Lowest);
  if (!expect_peek_and_advance(TokenTypes::RightParen)) {
    return nullptr;
  }
  return expr;
}

ast::ExprPtr Parser::parse_infix_expression(ast::ExprPtr left_expr) {
  const auto op = cur_token_.type;
  const auto prec = cur_precedence();
  next_token();

  auto right_expr = parse_expression(prec);

  return ast::make_expression<ast::InfixExpression>(std::move(left_expr), op,
                                                    std::move(right_expr));
}

std::vector<ast::ExprPtr>
Parser::parse_list_of_expressions(TokenTypes end_token) {
  // Should this consume the end-token or not?

  std::vector<ast::ExprPtr> expressions;
  if (peek_token_is(end_token)) {
    next_token();
    return expressions;
  }
  next_token();

  expressions.push_back(parse_expression(Precedence::Lowest));
  // should this be an input arg as well ...?
  while (peek_token_is(TokenTypes::Comma)) {
    next_token();
    next_token(); // Comma should be consumed
    expressions.push_back(parse_expression(Precedence::Lowest));
  }

  if (!expect_peek_and_advance(end_token)) {
    // may prefer a throw
    throw std::runtime_error("Unexpected Token at end of list " +
                             to_string(cur_token_));
  }
  return expressions;
}

ast::ExprPtr Parser::parse_function_expression(ast::ExprPtr func) {
  auto args = parse_list_of_expressions(TokenTypes::RightParen);
  return ast::make_expression<ast::FuncExpression>(std::move(func),
                                                   std::move(args));
}

ast::ExprPtr Parser::parse_array_expression() {
  return ast::make_expression<ast::ArrayExpression>(
      parse_list_of_expressions(TokenTypes::ArrayRightBracket));
}

Parser::Parser(Lexer lexer)
    : lexer_{std::move(lexer)},
      prefix_parse_fns_{{
          {TokenTypes::Identifier, &Parser::parse_identifier},
          {TokenTypes::Integer, &Parser::parse_integer_literal},
          {TokenTypes::Float, &Parser::parse_float_literal},
          {TokenTypes::String, &Parser::parse_string_literal},
          {TokenTypes::Minus, &Parser::parse_prefix_expression},
          {TokenTypes::Bang, &Parser::parse_prefix_expression},
          {TokenTypes::ArrayLeftBracket, &Parser::parse_array_expression},
      }},
      infix_parse_fns_{{
          {TokenTypes::Plus, &Parser::parse_infix_expression},
          {TokenTypes::Minus, &Parser::parse_infix_expression},
          {TokenTypes::Asterisk, &Parser::parse_infix_expression},
          {TokenTypes::Exp, &Parser::parse_infix_expression},
          {TokenTypes::Assign, &Parser::parse_infix_expression},
          {TokenTypes::Slash, &Parser::parse_infix_expression},
          {TokenTypes::Equal, &Parser::parse_infix_expression},
          {TokenTypes::NotEqual, &Parser::parse_infix_expression},
          {TokenTypes::GreaterThan, &Parser::parse_infix_expression},
          {TokenTypes::GreaterEqual, &Parser::parse_infix_expression},
          {TokenTypes::LessThan, &Parser::parse_infix_expression},
          {TokenTypes::LessEq, &Parser::parse_infix_expression},
          {TokenTypes::Or, &Parser::parse_infix_expression},
          {TokenTypes::And, &Parser::parse_infix_expression},
          {TokenTypes::Dot, &Parser::parse_infix_expression},
          {TokenTypes::LeftParen, &Parser::parse_function_expression},
      }} {
  next_token();
  next_token();
}

ast::ExprPtr Parser::parse() {
  // For now i'll assume we're parsing one expression statement at a time
  // and nothing more complicated
  auto expr = parse_expression(Precedence::Lowest);

  if (has_errors()) {
    throw ParserError(errors_);
  }
  return expr;
}

} // namespace edp::parser
