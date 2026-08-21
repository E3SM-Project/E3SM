#ifndef EDP_AST_HPP
#define EDP_AST_HPP

#include <edp/tokens.hpp>
#include <concepts>
#include <memory>
#include <string>
#include <variant>
#include <vector>

/**
 * @file ast.hpp
 * @brief Definition of AST Nodes, using std::variant
 */

namespace edp::ast {
/* A node of the AST is represented by an Expression which can be of any
 * type listed in the ExpressionVariant.
 * Expression contains a "visit" method that serves as a generic wrapper to
 * std::visit
 * */

struct Expression;

// Using a unique_ptr because i don't **think** there are multiple owners for
// one node
using ExprPtr = std::unique_ptr<const Expression>;

struct Identifier {
  std::string value;
};

struct PrefixExpression {
  TokenTypes op;
  ExprPtr right;
};

struct InfixExpression {
  ExprPtr left;
  TokenTypes op;
  ExprPtr right;
};

struct FuncExpression {
  ExprPtr function;
  std::vector<ExprPtr> args;
};

// struct BoundsExpression {
//   ExprPtr start;
//   ExprPtr stop;
// };

struct ArrayExpression {
  std::vector<ExprPtr> elements;
};

struct StringLiteral {
  std::string value;
};

struct FloatLiteral {
  float value;
};

struct IntegerLiteral {
  int value;
};

// Likely won't be needed.
// struct BooleanLiteral {
//   Token token;
//   bool value;
// };

using ExpressionVariant =
    std::variant<Identifier, PrefixExpression, InfixExpression, FuncExpression,
                 ArrayExpression, StringLiteral, FloatLiteral, IntegerLiteral>;

template <typename T>
concept ExpressionNode = std::constructible_from<ExpressionVariant, T&&>;

struct Expression {
  template <ExpressionNode T>
  explicit Expression(T&& value) : node_(std::forward<T>(value)) {}

  // This member function will be used for visitors so that
  // the node variant can remain private
  // NOTE: decltype(auto) allows visitors to return references
  template <typename Visitor> decltype(auto) visit(Visitor&& visitor) const {
    return std::visit(std::forward<Visitor>(visitor), node_);
  }

private:
  ExpressionVariant node_;
};

template <ExpressionNode Node, typename... Args>
  requires std::constructible_from<Node, Args&&...>
ExprPtr make_expression(Args&&... args) {
  return std::make_unique<const Expression>(Node{std::forward<Args>(args)...});
}

// Functions
std::string to_string(const Expression& expr);
// bool equal(const Expression& lhs, const Expression& rhs);

} // namespace edp::ast

#endif
