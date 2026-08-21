#include <charconv>
#include <dexpr/ast.hpp>
#include <dexpr/parser.hpp>
#include <dexpr/precedences.hpp>
#include <dexpr/tokens.hpp>
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
              to_string(peek_token_) + " at " + position_of(peek_token_));
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
    add_error("Illegal token " + to_string(peek_token_) + " at " +
              position_of(peek_token_));
  }
}

// prefix_parse_fns_{{
// }},
// infix_parse_fns_{{
// }}

Parser::PrefixFn Parser::get_prefix_parse_fn(TokenTypes tok_type) {
  switch (tok_type) {
  case TokenTypes::Identifier:
    return &Parser::parse_identifier;
  case TokenTypes::Integer:
    return &Parser::parse_integer_literal;
  case TokenTypes::Float:
    return &Parser::parse_float_literal;
  case TokenTypes::String:
    return &Parser::parse_string_literal;
  case TokenTypes::Minus:
  case TokenTypes::Bang:
    return &Parser::parse_prefix_expression;
  case TokenTypes::ArrayLeftBracket:
    return &Parser::parse_array_expression;
  case TokenTypes::LeftParen:
    return &Parser::parse_grouped_expression;
  default:
    return nullptr;
  }
}

Parser::InfixFn Parser::get_infix_parse_fn(TokenTypes tok_type) {
  switch (tok_type) {
  case TokenTypes::Plus:
  case TokenTypes::Minus:
  case TokenTypes::Asterisk:
  case TokenTypes::Exp:
  case TokenTypes::Assign:
  case TokenTypes::Slash:
  case TokenTypes::Equal:
  case TokenTypes::NotEqual:
  case TokenTypes::GreaterThan:
  case TokenTypes::GreaterEqual:
  case TokenTypes::LessThan:
  case TokenTypes::LessEq:
  case TokenTypes::Or:
  case TokenTypes::And:
  case TokenTypes::Dot:
    return &Parser::parse_infix_expression;
  case TokenTypes::LeftParen:
    return &Parser::parse_function_expression;
  default:
    return nullptr;
  }
}

/**
 * @brief Parses an expression using Pratt/operator-precedence parsing.
 *
 * Parsing begins by dispatching on the current token to construct the initial
 * left-hand expression using the appropriate prefix parse function. The parser
 * then examines subsequent tokens and extends that expression with infix
 * operations if the next operator has **higher** precedence.
 *
 * The Precedence prec argument sets the minimum binding precedence for
 * this invocation. Parsing stops when the next token has equal or lower
 * precedence, allowing the calling parse function to retain ownership of that
 * operator and thereby determine the grouping of the resulting AST.
 *
 * For example, when parsing:
 * @code
 * x + y * z
 * @endcode
 *
 * the multiplication operator has higher precedence than addition, so the
 * right-hand side of the addition is parsed as:
 * @code
 * y * z
 * @endcode
 *
 * producing the grouping:
 * @code
 * x + (y * z)
 * @endcode
 *
 * @param prec Minimum precedence required for subsequent operators to bind to
 *             the expression currently being parsed.
 * @return Pointer to the root expression node of the parsed AST subtree.
 * @throws ParserError If the current token cannot begin an expression.
 */
ast::ExprPtr Parser::parse_expression(Precedence prec) {
  const auto prefix_fn = get_prefix_parse_fn(cur_token_.type);
  if (prefix_fn == nullptr) {
    add_error("Unexpected Prefix Token " + to_string(cur_token_) + " at " +
              position_of(cur_token_));
    throw ParserError(errors_);
  }
  auto left_expr = (this->*prefix_fn)();

  while (!peek_token_is(TokenTypes::EndofFile) && prec < peek_precedence()) {
    const auto infix_fn = get_infix_parse_fn(peek_token_.type);
    if (infix_fn == nullptr) {
      return left_expr;
    }
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
  add_error("Integer literal out of range: " + cur_token_.literal + " at " +
            position_of(cur_token_));
  throw ParserError(errors_);
}

ast::ExprPtr Parser::parse_float_literal() {
  if (const auto value = parse_number<double>(cur_token_.literal)) {
    return ast::make_expression<ast::FloatLiteral>(*value);
  }
  add_error("Float literal out of range: " + cur_token_.literal + " at " +
            position_of(cur_token_));
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
  const auto prec = right_binding_precedence(op);
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

Parser::Parser(Lexer lexer) : lexer_{std::move(lexer)} {
  next_token();
  next_token();
}

/**
 * @brief Public entry point to parse an expression.
 *
 * Parses the string given to the Lexer used to construct the parser.
 *
 * @return Pointer to root of the AST.
 * @throws ParserError
 * @see Parser::parse_expression
 * @warning Assumes input has only **one** expression
 */
ast::ExprPtr Parser::parse() {
  auto expr = parse_expression(Precedence::Lowest);

  if (!peek_token_is(TokenTypes::EndofFile)) {
    add_error("Unexpected trailing input " + to_string(peek_token_) + " at " +
              position_of(peek_token_));
  }

  if (has_errors()) {
    throw ParserError(errors_);
  }
  return expr;
}

} // namespace dexpr::parser
