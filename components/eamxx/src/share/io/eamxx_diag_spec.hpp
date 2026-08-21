#ifndef SCREAM_EAMXX_DIAG_SPEC_HPP
#define SCREAM_EAMXX_DIAG_SPEC_HPP

#include <ekat_parameter_list.hpp>

#include <string>
#include <vector>

namespace scream {

/*
 * Lowering of a diagnostic expression into a tree of factory arguments.
 *
 * This is a *pure* step: it turns a string into data, and knows nothing about
 * grids, comms, field managers, or the diagnostic factory. That makes it
 * testable on its own, and keeps the "what does this name mean?" question in
 * one place instead of spread across a chain of regexes.
 *
 * Instantiating the tree (creating the diags, wiring inputs, computing an
 * evaluation order) is a separate step, and is not done here.
 */

/*
 * The names a node of the tree is known by.
 *
 *  - canonical: the fully-parenthesized print of the sub-expression this node
 *    represents. Used for two purely internal jobs: (a) it is the string a
 *    parent stores in its params (e.g. "field_name") to refer to this child,
 *    and which the child's field must therefore be aliased to, and (b) it is
 *    the key for reusing a diag that was already built for the same
 *    sub-expression. It must be deterministic and unambiguous. It is never
 *    shown to a user, and -- importantly -- never parsed back.
 *
 *  - registered: the name this node is known by in the FieldManager and, when
 *    written, in the output file. Empty for internal nodes. Because nothing
 *    ever parses it back, it does NOT need to be unambiguous: legacy names such
 *    as 'A_minus_B_over_C' are perfectly fine here even though they cannot be
 *    unambiguously re-read.
 *
 *  - write: whether the registered field is written to file. False for the
 *    intermediates declared in the 'aliases' section of the output yaml, which
 *    must exist in the FM but must not show up in the nc file.
 */
struct DiagNames {
  std::string canonical;
  std::string registered;
  bool        write = false;
};

struct DiagSpec {
  // Key to use in DiagnosticFactory. An EMPTY key marks an unresolved leaf:
  // names.canonical is a bare name, which is either a model field or a diag
  // registered under that very name. Only the FieldManager can tell the two
  // apart, so the choice is deliberately left to whoever builds the tree.
  std::string factory_key;

  // Arguments for the factory.
  // NOTE: "grid_name" is deliberately NOT set here, since lowering is
  //       grid-independent. Whoever instantiates the tree must inject it.
  //       Diags that do not use it simply ignore the extra parameter.
  ekat::ParameterList params;

  // The operands of this node, in the order they appear in the expression.
  // A node may also need inputs that never appear in the expression (e.g.
  // VertContract with dp weighting needs 'pseudo_density'); those are not
  // children, and must be resolved against the FM at instantiation time.
  std::vector<DiagSpec> children;

  DiagNames names;

  bool is_leaf () const { return factory_key.empty(); }
};

// Parse 'expr' with the e3sm diags parser, and lower the resulting AST into a
// DiagSpec tree. Throws if 'expr' is not a valid diagnostic expression.
// 'registered'/'write' are attached to the root node (see DiagNames above);
// all other nodes are internal, and get a canonical name only.
DiagSpec lower_to_diag_spec (const std::string& expr,
                             const std::string& registered = "",
                             const bool write = true);

} // namespace scream

#endif // SCREAM_EAMXX_DIAG_SPEC_HPP
