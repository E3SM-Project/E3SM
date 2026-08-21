#include <edp/ast.hpp>
#include <edp/tokens.hpp>
#include <algorithm>
#include <cmath>
#include <format>
#include <span>

/**
 * @file ast_print.cpp
 * @brief Implementation of Visitors for printing AST nodes
 */
namespace edp::ast {

namespace {

struct ToStringVisitor {

  std::string operator()(const Identifier& expr) const;
  std::string operator()(const PrefixExpression& expr) const;
  std::string operator()(const InfixExpression& expr) const;
  std::string operator()(const FuncExpression& expr) const;
  std::string operator()(const ArrayExpression& expr) const;
  std::string operator()(const StringLiteral& expr) const;
  std::string operator()(const FloatLiteral& expr) const;
  std::string operator()(const IntegerLiteral& expr) const;
};

std::string expr_list_to_string(std::span<const ExprPtr> vals) {

  std::string result;
  bool first = true;

  std::ranges::for_each(vals, [&](const ExprPtr& val) {
    if (!first) {
      result += ", ";
    }
    first = false;
    result += to_string(*val);
  });
  return result;
}

std::string ToStringVisitor::operator()(const Identifier& expr) const {
  return expr.value;
};

std::string ToStringVisitor::operator()(const PrefixExpression& expr) const {
  return "(" + unary_op_to_string(expr.op) + to_string(*expr.right) + ")";
};

std::string ToStringVisitor::operator()(const InfixExpression& expr) const {
  // '.' binds at Call, the tightest level, so wrapping it adds noise without
  // disambiguating anything...
  if (expr.op == TokenTypes::Dot) {
    return to_string(*expr.left) + "." + to_string(*expr.right);
  }
  return "(" + to_string(*expr.left) + binary_op_to_string(expr.op) +
         to_string(*expr.right) + ")";
};
std::string ToStringVisitor::operator()(const FuncExpression& expr) const {
  return to_string(*expr.function) + "(" + expr_list_to_string(expr.args) + ")";
};

std::string ToStringVisitor::operator()(const ArrayExpression& expr) const {
  return +"[" + expr_list_to_string(expr.elements) + "]";
};

std::string ToStringVisitor::operator()(const StringLiteral& expr) const {
  return "'" + expr.value + "'";
};
std::string ToStringVisitor::operator()(const IntegerLiteral& expr) const {
  return std::to_string(expr.value);
};
std::string ToStringVisitor::operator()(const FloatLiteral& expr) const {
  auto result = std::format("{}", expr.value);
  if (std::isfinite(expr.value) &&
      result.find_first_of(".eE") == std::string::npos) {
    result += ".0";
  }
  return result;
};

} // namespace

std::string to_string(const Expression& expr) {
  return expr.visit(ToStringVisitor{});
}

} // namespace edp::ast
