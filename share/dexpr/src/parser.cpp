#include <dexpr/parser.hpp>
#include <dexpr/ast.hpp>
#include <dexpr/precedences.hpp>
#include <dexpr/tokens.hpp>
#include <charconv>
#include <optional>
#include <stdexcept>
#include <string>

namespace dexpr::parser {

namespace {

// std::from_chars rather than std::sto*: it is locale-independent (std::stof
// reads "1.5" as 1 wherever ',' is the decimal separator), it does not throw,
// and it reports where it stopped. Requiring it to stop at the end of the
// literal means a malformed one is rejected outright instead of being read as
// its leading prefix, which is exactly how "0.5.3" used to become 0.5.
template <typename T>
std::optional<T> parse_number(const std::string& literal) {
  T value{};
  const auto* const first = literal.data();
  const auto* const last = first + literal.size();

  const auto [stopped_at, ec] = std::from_chars(first, last, value);
  if (ec != std::errc{} || stopped_at != last) {
    return std::nullopt;
  }
  return value;
}

} // namespace

bool Parser::peek_token_is(TokenTypes expected_type) {
  return peek_token_.type == expected_type;
};

void Parser::add_error(std::string msg) { errors_.push_back(std::move(msg)); }

bool Parser::expect_peek_and_advance(TokenTypes expected_type) {
  if (peek_token_is(expected_type)) {
    next_token();
    return true;
  } else {
    add_error("Expected " + std::string(to_string(expected_type)) + ", got " +
              to_string(peek_token_));
    return false;
  }
}

Precedence Parser::peek_precedence() {
  return token_precedence(peek_token_.type);
}

bool Parser::has_errors() { return !errors_.empty(); }

void Parser::next_token() {
  cur_token_ = peek_token_;
  peek_token_ = lexer_.next_token();
  if (peek_token_is(TokenTypes::Illegal)) {
    add_error("Illegal token " + to_string(peek_token_));
  }
}

ast::ExprPtr Parser::parse_expression(Precedence prec) {
  const auto prefix = prefix_parse_fns_.find(cur_token_.type);
  if (prefix == prefix_parse_fns_.end()) {
    add_error("Unexpected Prefix Token " + to_string(cur_token_));
    throw ParserError(errors_);
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
  if (const auto value = parse_number<int>(cur_token_.literal)) {
    return ast::make_expression<ast::IntegerLiteral>(*value);
  }
  add_error("Integer literal out of range: " + cur_token_.literal);
  throw ParserError(errors_);
}

ast::ExprPtr Parser::parse_float_literal() {
  if (const auto value = parse_number<double>(cur_token_.literal)) {
    return ast::make_expression<ast::FloatLiteral>(*value);
  }
  add_error("Float literal out of range: " + cur_token_.literal);
  throw ParserError(errors_);
}
ast::ExprPtr Parser::parse_prefix_expression() {
  auto op = cur_token_.type;
  next_token();
  auto right_expr = parse_expression(Precedence::Prefix);
  return ast::make_expression<ast::PrefixExpression>(op, std::move(right_expr));
}

ast::ExprPtr Parser::parse_grouped_expression() {
  next_token();
  auto expr = parse_expression(Precedence::Lowest);
  if (!expect_peek_and_advance(TokenTypes::RightParen)) {
    // Throwing rather than returning null: parse() would catch this at the end
    // anyway, but until then the null travels through the tree builders, and
    // every visitor dereferences its children unguarded.
    throw ParserError(errors_);
  }
  return expr;
}

ast::ExprPtr Parser::parse_infix_expression(ast::ExprPtr left_expr) {
  const auto op = cur_token_.type;
  const auto prec = cur_precedence(op);
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
    throw ParserError(errors_);
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
          {TokenTypes::LeftParen, &Parser::parse_grouped_expression},
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

  // TODO: maybe this is a bad idea?
  // The expression must account for the whole input
  if (!peek_token_is(TokenTypes::EndofFile)) {
    add_error("Unexpected trailing input at " + to_string(peek_token_));
  }

  if (has_errors()) {
    throw ParserError(errors_);
  }
  return expr;
}

} // namespace dexpr::parser
