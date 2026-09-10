/**
 * @file precedences.hpp
 * @brief Defines operator precedence for diagnostic expressions.
 *
 * Precedence determines how tightly operators bind to their operands during
 * expression parsing. Higher precedence values bind more tightly than lower
 * precedence values.
 *
 * Each operator TokenType is mapped to a Precedence level by
 * token_precedence(). The Pratt parser uses these levels to determine the
 * grouping of expressions without requiring explicit parentheses.
 *
 * For example, given:
 * @code
 * x = -w * y
 * @endcode
 *
 * with:
 * @code
 * Precedence::Equal < Precedence::Product < Precedence::Prefix
 * @endcode
 *
 * the expression is parsed as:
 * @code
 * x = ((-w) * y)
 * @endcode
 */
#ifndef DEXPR_PRECEDENCES_HPP
#define DEXPR_PRECEDENCES_HPP

#include <dexpr/tokens.hpp>
namespace dexpr::parser {

enum class Precedence {
  Lowest,
  Assignment,
  Logical,
  Equalitative,
  Comparison,
  Additive,
  Multiplicative,
  Prefix,
  Exponent,
  Call,
};

Precedence token_precedence(TokenTypes type);
Precedence right_binding_precedence(TokenTypes type);
} // namespace dexpr::parser

#endif
