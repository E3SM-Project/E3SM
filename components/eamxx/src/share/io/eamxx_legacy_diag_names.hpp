#ifndef SCREAM_EAMXX_LEGACY_DIAG_NAMES_HPP
#define SCREAM_EAMXX_LEGACY_DIAG_NAMES_HPP

#include <string>

namespace scream {

/*
 * Translation of the legacy underscore-mangled diagnostic names into the
 * expression syntax.
 *
 * This is a compatibility layer, and is meant to be deleted once the yaml files
 * in the wild have moved to expressions. Everything it knows lives in this one
 * file, so that removing it is a single commit.
 *
 * Two properties matter:
 *
 *  - it is NOT recursive. Only the outermost pattern is translated, and the
 *    operands are left as bare legacy names: 'A_minus_B_over_C' becomes
 *    'A_minus_B/C', where 'A_minus_B' is still a single identifier. This is
 *    deliberate. The output layer has always checked the FieldManager before
 *    treating a name as a diagnostic, so a model field genuinely called
 *    'foo_prev' must stay a field; translating operands eagerly would rewrite
 *    it to 'foo.prev()' and silently compute something else. Recursion happens
 *    where the FieldManager is in reach, i.e. while walking the tree.
 *
 *  - the order the patterns are tried in mirrors, exactly, the order of the
 *    regexes it replaces. That order is load-bearing: '_over_dt' must be tried
 *    before the binary ops so that 'X_over_dt' is not read as X over dt, and
 *    the binary ops must be tried before '_prev' so that 'X_minus_X_prev' is
 *    a subtraction rather than the previous value of 'X_minus_X'.
 *
 * Returns the input unchanged if it matches no legacy pattern, which is how the
 * caller can tell a translation happened.
 */
std::string legacy_to_expr (const std::string& name);

} // namespace scream

#endif // SCREAM_EAMXX_LEGACY_DIAG_NAMES_HPP
